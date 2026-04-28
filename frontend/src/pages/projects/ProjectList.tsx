/**
 * Project List Page — TanStack Query migration.
 */
import type React from 'react'
import { useMemo, useState } from 'react'
import { Button, Space, Dropdown, Modal } from 'antd'
import {
  PlusOutlined,
  FolderOutlined,
  EditOutlined,
  DeleteOutlined,
  MoreOutlined,
  InboxOutlined,
  RestOutlined,
  CheckCircleOutlined,
  ClockCircleOutlined,
} from '@ant-design/icons'
import type { ColumnsType } from 'antd/es/table'
import dayjs from 'dayjs'
import relativeTime from 'dayjs/plugin/relativeTime'

import {
  useProjectList,
  useGlobalProjectStats,
  useDeleteProject,
  useArchiveProject,
  useRestoreProject,
} from '@/hooks/queries'
import { ProjectFormModal } from './components/ProjectFormModal'
import {
  PageHeader,
  DataTable,
  StatisticCard,
  StatusTag,
  FilterBar,
  EnhancedEmptyState,
  PageSkeleton,
  FadeIn,
  StaggeredList,
} from '@/components/common'
import type { FilterConfig, StatisticItem } from '@/components/common'
import { toast, notifications } from '@/utils/notification'
import { useFilters } from '@/hooks'
import type { Project } from '@/types/project'

dayjs.extend(relativeTime)

export const ProjectList: React.FC = () => {
  const { data: listData, isLoading, isFetching } = useProjectList()
  const { data: stats } = useGlobalProjectStats()

  const items: Project[] = useMemo(() => (listData as any)?.items ?? (listData as any) ?? [], [listData])

  const [modalOpen, setModalOpen] = useState(false)
  const [editingProject, setEditingProject] = useState<Project | null>(null)

  const deleteMutation = useDeleteProject()
  const archiveMutation = useArchiveProject()
  const restoreMutation = useRestoreProject()

  const {
    filters,
    setFilter,
    resetFilters: handleFilterReset,
  } = useFilters({
    initialFilters: { search: '', status: 'all' },
  })

  const filterConfigs: FilterConfig[] = [
    { type: 'search', key: 'search', placeholder: 'Search items...' },
    {
      type: 'select',
      key: 'status',
      label: 'Status',
      options: [
        { label: 'All Status', value: 'all' },
        { label: 'Active', value: 'active' },
        { label: 'Archived', value: 'archived' },
        { label: 'Completed', value: 'completed' },
      ],
    },
  ]

  const filteredProjects = useMemo(
    () =>
      items.filter((project) => {
        const matchesStatus = filters.status === 'all' || project.status === filters.status
        const matchesSearch =
          filters.search === '' ||
          project.name.toLowerCase().includes(filters.search.toLowerCase()) ||
          project.description?.toLowerCase().includes(filters.search.toLowerCase())
        return matchesStatus && matchesSearch
      }),
    [items, filters],
  )

  const handleCreate = () => {
    setEditingProject(null)
    setModalOpen(true)
  }

  const handleEdit = (project: Project) => {
    setEditingProject(project)
    setModalOpen(true)
  }

  const handleDelete = (project: Project) => {
    Modal.confirm({
      title: '删除项目',
      content: `您确定要删除项目 "${project.name}" 吗？这将同时删除所有关联的样本、文件和任务。`,
      okText: '确认',
      cancelText: '取消',
      okType: 'danger',
      onOk: async () => {
        const loadingToast = toast.loading('删除中...')
        try {
          await deleteMutation.mutateAsync(project.id)
          loadingToast()
          notifications.deleteSuccess()
        } catch {
          loadingToast()
          notifications.deleteError()
        }
      },
    })
  }

  const handleArchive = async (project: Project) => {
    const loadingToast = toast.loading('归档中...')
    try {
      await archiveMutation.mutateAsync(project.id)
      loadingToast()
      toast.success('项目已归档')
    } catch {
      loadingToast()
      toast.error('归档失败，请重试')
    }
  }

  const handleRestore = async (project: Project) => {
    const loadingToast = toast.loading('恢复中...')
    try {
      await restoreMutation.mutateAsync(project.id)
      loadingToast()
      toast.success('项目已恢复')
    } catch {
      loadingToast()
      toast.error('恢复失败，请重试')
    }
  }

  const columns: ColumnsType<Project> = [
    {
      title: 'Project Name',
      dataIndex: 'name',
      key: 'name',
      render: (name, record) => (
        <Space>
          <FolderOutlined style={{ fontSize: 18, color: 'var(--color-primary)' }} />
          <span style={{ fontWeight: 500 }}>{name}</span>
          {record.status === 'archived' && <StatusTag status="archived" />}
        </Space>
      ),
    },
    { title: 'Description', dataIndex: 'description', key: 'description', ellipsis: true },
    {
      title: 'Status',
      dataIndex: 'status',
      key: 'status',
      width: 120,
      render: (status) => <StatusTag status={status} />,
    },
    {
      title: 'Samples',
      dataIndex: 'sample_count',
      key: 'sample_count',
      width: 100,
      render: (count) => count || 0,
    },
    {
      title: 'Tasks',
      dataIndex: 'task_count',
      key: 'task_count',
      width: 100,
      render: (count) => count || 0,
    },
    {
      title: 'Created',
      dataIndex: 'created_at',
      key: 'created_at',
      width: 150,
      render: (date) => dayjs(date).fromNow(),
    },
    {
      title: 'Actions',
      key: 'actions',
      width: 100,
      render: (_, record) => (
        <Dropdown
          menu={{
            items: [
              { key: 'edit', label: 'Edit', icon: <EditOutlined />, onClick: () => handleEdit(record) },
              {
                key: 'archive',
                label: record.status === 'archived' ? 'Restore' : 'Archive',
                icon: record.status === 'archived' ? <RestOutlined /> : <InboxOutlined />,
                onClick: () => (record.status === 'archived' ? handleRestore(record) : handleArchive(record)),
              },
              { type: 'divider' },
              {
                key: 'delete',
                label: 'Delete',
                icon: <DeleteOutlined />,
                danger: true,
                onClick: () => handleDelete(record),
              },
            ],
          }}
        >
          <Button icon={<MoreOutlined />} />
        </Dropdown>
      ),
    },
  ]

  const statisticItems: StatisticItem[] = [
    {
      key: 'total',
      title: 'Total Projects',
      value: (stats as any)?.total_projects || 0,
      prefix: <FolderOutlined />,
      valueStyle: { color: 'var(--color-primary)' },
    },
    {
      key: 'active',
      title: 'Active Projects',
      value: (stats as any)?.active_projects || 0,
      prefix: <CheckCircleOutlined />,
      valueStyle: { color: 'var(--color-success)' },
    },
    {
      key: 'total_tasks',
      title: 'Total Tasks',
      value: (stats as any)?.total_tasks || 0,
      prefix: <ClockCircleOutlined />,
      valueStyle: { color: '#722ed1' },
    },
    {
      key: 'active_tasks',
      title: 'Active Tasks',
      value: (stats as any)?.active_tasks || 0,
      prefix: <ClockCircleOutlined />,
      valueStyle: { color: 'var(--color-warning)' },
    },
  ]

  if (isLoading && items.length === 0) {
    return <PageSkeleton hasHeader hasFilters rows={8} />
  }

  return (
    <div>
      <StaggeredList staggerDelay={80} baseDelay={0} direction="up">
        <StatisticCard items={statisticItems} />
      </StaggeredList>

      <FadeIn direction="up" delay={100} duration={300}>
        <PageHeader
          left={
            <FilterBar
              filters={filterConfigs}
              values={filters}
              onFilterChange={(key, value) => setFilter(key as 'search' | 'status', value)}
              onReset={handleFilterReset}
            />
          }
          right={
            <Button type="primary" icon={<PlusOutlined />} onClick={handleCreate}>
              New Project
            </Button>
          }
        />
      </FadeIn>

      <FadeIn direction="up" delay={200} duration={300}>
        {filteredProjects.length === 0 && !isFetching ? (
          <EnhancedEmptyState
            type={filters.search || filters.status !== 'all' ? 'noSearchResults' : 'noData'}
            title={filters.search || filters.status !== 'all' ? 'No matching items' : 'No items yet'}
            description={
              filters.search || filters.status !== 'all'
                ? 'Try adjusting your search criteria or filters'
                : 'Create your first project to get started with NGS analysis'
            }
            action={{ text: 'Create Project', onClick: handleCreate, icon: <PlusOutlined /> }}
            size="default"
          />
        ) : (
          <DataTable
            columns={columns}
            dataSource={filteredProjects}
            rowKey="id"
            loading={isFetching && items.length === 0}
            pagination={{
              pageSize: 10,
              showSizeChanger: true,
              showTotal: (total) => `Total ${total} items`,
            }}
            emptyText="No Projects"
            emptyDescription="Create your first project to get started"
          />
        )}
      </FadeIn>

      <ProjectFormModal
        open={modalOpen}
        project={editingProject}
        onClose={() => {
          setModalOpen(false)
          setEditingProject(null)
        }}
        onSuccess={() => {
          setModalOpen(false)
          setEditingProject(null)
          // mutations inside the modal already invalidate the cache
        }}
      />
    </div>
  )
}

export default ProjectList
