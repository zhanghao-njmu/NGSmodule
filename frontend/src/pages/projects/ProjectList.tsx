/**
 * Project List Page - Complete project management interface
 */
import type React from 'react'
import { useEffect, useState } from 'react'
import { Button, Space, Dropdown } from 'antd'
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
import { useProjectStore } from '../../store/projectStore'
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
} from '../../components/common'
import type { FilterConfig } from '../../components/common'
import { toast, notifications } from '../../utils/notification'
import { useFilters } from '@/hooks'
import type { StatisticItem } from '../../components/common'
import type { Project } from '../../types/project'
import dayjs from 'dayjs'
import relativeTime from 'dayjs/plugin/relativeTime'

dayjs.extend(relativeTime)

export const ProjectList: React.FC = () => {
  const { items, stats, loading, fetchItems, fetchStats, deleteItem, archiveProject, restoreProject } =
    useProjectStore()

  const [modalOpen, setModalOpen] = useState(false)
  const [editingProject, setEditingProject] = useState<Project | null>(null)

  // Using useFilters hook eliminates repetitive filter state management
  const {
    filters,
    setFilter,
    resetFilters: handleFilterReset,
  } = useFilters({
    initialFilters: {
      search: '',
      status: 'all',
    },
  })

  useEffect(() => {
    fetchItems()
    fetchStats()
  }, [])

  // Filter configuration for FilterBar
  const filterConfigs: FilterConfig[] = [
    {
      type: 'search',
      key: 'search',
      placeholder: 'Search items...',
    },
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

  // Filter items based on status and search
  const filteredProjects = items.filter((project) => {
    const matchesStatus = filters.status === 'all' || project.status === filters.status
    const matchesSearch =
      filters.search === '' ||
      project.name.toLowerCase().includes(filters.search.toLowerCase()) ||
      project.description?.toLowerCase().includes(filters.search.toLowerCase())
    return matchesStatus && matchesSearch
  })

  const handleCreate = () => {
    setEditingProject(null)
    setModalOpen(true)
  }

  const handleEdit = (project: Project) => {
    setEditingProject(project)
    setModalOpen(true)
  }

  const handleDelete = (project: Project) => {
    // Use custom dangerous action confirmation for critical operations
    confirmDangerousAction(
      '删除项目',
      `您确定要删除项目 "${project.name}" 吗？这将同时删除所有关联的样本、文件和任务。`,
      async () => {
        const loadingToast = toast.loading('删除中...')
        try {
          await deleteItem(project.id)
          loadingToast()
          notifications.deleteSuccess()
          fetchItems() // Refresh list
          fetchStats() // Refresh stats
        } catch (error) {
          loadingToast()
          notifications.deleteError()
        }
      },
    )
  }

  const handleArchive = async (project: Project) => {
    const loadingToast = toast.loading('归档中...')
    try {
      await archiveProject(project.id)
      loadingToast()
      toast.success('项目已归档')
      fetchItems()
      fetchStats()
    } catch (error) {
      loadingToast()
      toast.error('归档失败，请重试')
    }
  }

  const handleRestore = async (project: Project) => {
    const loadingToast = toast.loading('恢复中...')
    try {
      await restoreProject(project.id)
      loadingToast()
      toast.success('项目已恢复')
      fetchItems()
      fetchStats()
    } catch (error) {
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
    {
      title: 'Description',
      dataIndex: 'description',
      key: 'description',
      ellipsis: true,
    },
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
              {
                key: 'edit',
                label: 'Edit',
                icon: <EditOutlined />,
                onClick: () => handleEdit(record),
              },
              {
                key: 'archive',
                label: record.status === 'archived' ? 'Restore' : 'Archive',
                icon: record.status === 'archived' ? <RestOutlined /> : <InboxOutlined />,
                onClick: () => (record.status === 'archived' ? handleRestore(record) : handleArchive(record)),
              },
              {
                type: 'divider',
              },
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
      value: stats?.total_projects || 0,
      prefix: <FolderOutlined />,
      valueStyle: { color: 'var(--color-primary)' },
    },
    {
      key: 'active',
      title: 'Active Projects',
      value: stats?.active_projects || 0,
      prefix: <CheckCircleOutlined />,
      valueStyle: { color: 'var(--color-success)' },
    },
    {
      key: 'total_tasks',
      title: 'Total Tasks',
      value: stats?.total_tasks || 0,
      prefix: <ClockCircleOutlined />,
      valueStyle: { color: '#722ed1' },
    },
    {
      key: 'active_tasks',
      title: 'Active Tasks',
      value: stats?.active_tasks || 0,
      prefix: <ClockCircleOutlined />,
      valueStyle: { color: 'var(--color-warning)' },
    },
  ]

  // Show skeleton while loading
  if (loading && items.length === 0) {
    return <PageSkeleton hasHeader hasFilters rows={8} />
  }

  return (
    <div>
      {/* Statistics Cards with Fade In Animation */}
      <FadeIn direction="up" delay={0} duration={300}>
        <StatisticCard items={statisticItems} />
      </FadeIn>

      {/* Filters and Actions */}
      <FadeIn direction="up" delay={100} duration={300}>
        <PageHeader
          left={
            <FilterBar
              filters={filterConfigs}
              values={filters}
              onFilterChange={setFilter}
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

      {/* Projects Table with Fade In Animation */}
      <FadeIn direction="up" delay={200} duration={300}>
        {filteredProjects.length === 0 && !loading ? (
          <EnhancedEmptyState
            type={filters.search || filters.status !== 'all' ? 'noSearchResults' : 'noData'}
            title={filters.search || filters.status !== 'all' ? 'No matching items' : 'No items yet'}
            description={
              filters.search || filters.status !== 'all'
                ? 'Try adjusting your search criteria or filters'
                : 'Create your first project to get started with NGS analysis'
            }
            action={{
              text: 'Create Project',
              onClick: handleCreate,
              icon: <PlusOutlined />,
            }}
            size="default"
          />
        ) : (
          <DataTable
            columns={columns}
            dataSource={filteredProjects}
            rowKey="id"
            loading={loading}
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

      {/* Create/Edit Modal */}
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
          fetchItems()
          fetchStats()
        }}
      />
    </div>
  )
}

export default ProjectList
