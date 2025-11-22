/**
 * Project List Page - Complete project management interface
 */
import React, { useEffect, useState } from 'react'
import {
  Button,
  Space,
  Dropdown,
  Input,
  Select,
} from 'antd'
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
import { PageHeader, DataTable, StatisticCard, StatusTag } from '../../components/common'
import { confirmDelete, confirmDangerousAction } from '../../components/common/ConfirmDialog'
import { toast, notifications } from '../../utils/notification'
import type { StatisticItem } from '../../components/common'
import type { Project } from '../../types/project'
import dayjs from 'dayjs'
import relativeTime from 'dayjs/plugin/relativeTime'

dayjs.extend(relativeTime)

const { Search } = Input
const { Option } = Select

export const ProjectList: React.FC = () => {
  const {
    projects,
    stats,
    loading,
    fetchProjects,
    fetchStats,
    deleteProject,
    archiveProject,
    restoreProject,
  } = useProjectStore()

  const [modalOpen, setModalOpen] = useState(false)
  const [editingProject, setEditingProject] = useState<Project | null>(null)
  const [statusFilter, setStatusFilter] = useState<string>('all')
  const [searchText, setSearchText] = useState('')

  useEffect(() => {
    fetchProjects()
    fetchStats()
  }, [])

  // Filter projects based on status and search
  const filteredProjects = projects.filter((project) => {
    const matchesStatus = statusFilter === 'all' || project.status === statusFilter
    const matchesSearch =
      searchText === '' ||
      project.name.toLowerCase().includes(searchText.toLowerCase()) ||
      project.description?.toLowerCase().includes(searchText.toLowerCase())
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
          await deleteProject(project.id)
          loadingToast()
          notifications.deleteSuccess()
          fetchProjects() // Refresh list
          fetchStats() // Refresh stats
        } catch (error) {
          loadingToast()
          notifications.deleteError()
        }
      }
    )
  }

  const handleArchive = async (project: Project) => {
    const loadingToast = toast.loading('归档中...')
    try {
      await archiveProject(project.id)
      loadingToast()
      toast.success('项目已归档')
      fetchProjects()
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
      fetchProjects()
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
                onClick: () =>
                  record.status === 'archived'
                    ? handleRestore(record)
                    : handleArchive(record),
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

  return (
    <div>
      {/* Statistics Cards */}
      <StatisticCard items={statisticItems} />

      {/* Filters and Actions */}
      <PageHeader
        left={
          <>
            <Search
              placeholder="Search projects..."
              style={{ width: 300 }}
              value={searchText}
              onChange={(e) => setSearchText(e.target.value)}
              allowClear
            />
            <Select
              value={statusFilter}
              onChange={setStatusFilter}
              style={{ width: 150 }}
            >
              <Option value="all">All Status</Option>
              <Option value="active">Active</Option>
              <Option value="archived">Archived</Option>
              <Option value="completed">Completed</Option>
            </Select>
          </>
        }
        right={
          <Button type="primary" icon={<PlusOutlined />} onClick={handleCreate}>
            New Project
          </Button>
        }
      />

      {/* Projects Table */}
      <DataTable
        columns={columns}
        dataSource={filteredProjects}
        rowKey="id"
        loading={loading}
        pagination={{
          pageSize: 10,
          showSizeChanger: true,
          showTotal: (total) => `Total ${total} projects`,
        }}
        emptyText="No Projects"
        emptyDescription="Create your first project to get started"
      />

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
          fetchProjects()
          fetchStats()
        }}
      />
    </div>
  )
}

export default ProjectList
